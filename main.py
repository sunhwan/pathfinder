from flask import Flask, render_template, request, redirect, url_for, session, abort, send_file
from werkzeug import secure_filename
from flask import jsonify
import os
import filecmp
import uuid
import sqlite3
from flask import g
from datetime import datetime
import conf
import subprocess as sp

ALLOWED_EXTENSIONS = set(['pdb'])

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = conf.UPLOAD_FOLDER
app.config['JOB_FOLDER'] = conf.JOB_FOLDER

def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        is_new = not os.path.exists(conf.DATABASE)
        db = g._database = sqlite3.connect(conf.DATABASE)
        # initialize?
        if is_new: init_db()
    return db 

def init_db():
    with app.app_context():
        db = get_db()
        with app.open_resource('schema.sql', mode='r') as f:
            db.cursor().executescript(f.read())
        db.commit()

def get_job_folder(uuid):
    return os.path.join(conf.JOB_FOLDER, uuid)

def check_taskserver_status():
    if conf.TASKSERVER == 'pbs':
        try:        
            sp.check_output('qstat')
            return True
        except:
            return False
    else:
        import socket
        try:
            ip, port = open(os.path.join(conf.JOB_FOLDER, 'taskserver.info')).readline().strip().split(':')
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.connect((ip, int(port)))
            sock.sendall('ok:')
            response = sock.recv(1024)
            if response == 'OK': return True
        except:
            pass

        print 'create taskserver...'
        sp.Popen(['python', 'taskserver.py'], cwd=conf.BASEDIR)

    return True

def submit_task(jobid, email):
    db = get_db()
    cur = db.cursor()
    jobdir = os.path.join(app.config['JOB_FOLDER'], jobid)
    cur.execute('insert into job( uuid, email, date, status ) values (?, ?, ?, ?)', (jobid, email, datetime.now(), 'Q'))
    db.commit()
    
    os.environ['USER'] = 'sunhwan'
    p = sp.Popen(['/usr/torque/bin/qsub', 'run.pbs'], cwd=jobdir, stdout=sp.PIPE, stderr=sp.PIPE)
    pid = p.communicate()[0]

    from flask_mail import Mail, Message
    mail = Mail(app)
    msg = Message('[ANMPathway] Your job is submitted.', sender=conf.ADMIN_EMAIL, recipients=[email], bcc=[conf.ADMIN_EMAIL, conf.EXTRA_EMAIL])
    msg.html = render_template('email_submit.tpl', uuid=jobid)
    mail.send(msg)

    return True


# main

@app.teardown_appcontext
def close_connection(exception):
    db = getattr(g, '_database', None)
    if db is not None:
        db.close()

@app.errorhandler(404)
def page_not_found(e):
    return render_template('404.html'), 404

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/tutorial')
def tutorial():
    return render_template('tutorial.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/download/<uuid>/<filename>')
@app.route('/download/<uuid>', defaults={'filename': None})
def download(uuid, filename):
    jobdir = get_job_folder(uuid)
    if not os.path.exists(jobdir): abort(404)

    if filename:
        filename = os.path.join(jobdir, os.path.basename(filename))
        if not os.path.exists(filename): abort(404)
        fp = open(filename)
        if filename.endswith('gif'):
            return send_file(fp, as_attachment=False, attachment_filename=os.path.basename(filename))
        return send_file(fp, as_attachment=True, attachment_filename=os.path.basename(filename))

    else:
        import tarfile
        import StringIO, glob
        fp = StringIO.StringIO()
        tar = tarfile.open(fileobj=fp, mode='w:gz')
        excludes = ['taskmanager.py', 'err', 'out', 'run.pbs']
        for f in glob.glob('%s/*' % jobdir):
            if os.path.basename(f) in excludes: continue
            tar.add(f, arcname='anmpathway/%s' % os.path.basename(f))
        tar.close()
        fp.seek(0)
        return send_file(fp, mimetype='application/x-gzip', as_attachment=True, attachment_filename='anmpathway.tar.gz')

@app.route('/status/<uuid>')
def status(uuid):
    cur = get_db().cursor()
    cur.execute('select * from job where uuid=?', (uuid,))
    row = cur.fetchone()
    uuid, email, date, status = row

    check_taskserver_status()
    return render_template('queue.html', uuid=uuid, date=date, status=status)

@app.route('/success/<uuid>')
def success(uuid):
    return render_template('success.html', uuid=uuid)

@app.route('/', methods=['GET', 'POST'])
def pathfinder():
    if request.method == 'GET': return redirect(url_for('index'))

    # sanity check
    flag = False
    message = {590: 'Missing PDB file(s)',
               591: 'Please select at least one chain',
               592: 'Please provide e-mail address', 
               593: 'Please provide proper residue range',
               594: 'Some parameters are missing for initial and final PDB'}

    if not os.path.exists(os.path.join(app.config['UPLOAD_FOLDER'], request.form.get('pdb1', ''))): flag = 590
    if not os.path.exists(os.path.join(app.config['UPLOAD_FOLDER'], request.form.get('pdb2', ''))): flag = 590
    if not request.form.getlist('chains'): flag = 591
    if not request.form.getlist('email'): flag = 592

    for chain in request.form.getlist('chains'):
        resid1 = map(int, request.form.getlist('chain.%s.from' % chain))
        resid2 = map(int, request.form.getlist('chain.%s.to' % chain))
        if len(resid1) != len(resid2): flag = 593

    for pname in ['fc1', 'fc2', 'cutoff1', 'cutoff2', 'offset1', 'offset2']:
        if not request.form.get(pname, None): flag = 594

    if flag:
        #error = {'code': 599, 'message': 'Something went terribly wrong'}
        error = {'code': flag, 'message': message[flag]}
        return render_template('index.html', error=error)

    jobid = str(uuid.uuid4())
    jobdir = os.path.join(app.config['JOB_FOLDER'], jobid)
    os.mkdir(jobdir)

    pdb1 = secure_filename(request.form['pdb1'])
    pdb2 = secure_filename(request.form['pdb2'])
    pdb1_uploaded = os.path.join(app.config['UPLOAD_FOLDER'], pdb1)
    pdb2_uploaded = os.path.join(app.config['UPLOAD_FOLDER'], pdb2)
    email = request.form['email']
    chains = request.form.getlist('chains')

    from shutil import copy
    copy(pdb1_uploaded, os.path.join(jobdir, pdb1))
    copy(pdb2_uploaded, os.path.join(jobdir, pdb2))

    executables = ['static/src/step-1-find-initial-point/1_locate_struct_on_cusp_v2', 
                   'static/src/step-2-minimize-on-cusp/2_find_min_on_cusp_v4', 
                   'static/src/step-3-slide-down/3_desc_one_surface_v3',
                   'static/src/step-4-collect-structures/4_collec_ener_v2',
                   'static/src/step-5-make-pathway/5_makepathway_v2',
                   'static/src/prepare_input_structure_files.tcl',
                   'static/src/create_input_files_ANMPathway.pl',
                   'static/src/non-native-contacts/find_pairs_path_v2',
                   'static/src/taskmanager.py',
                   'static/src/movie.tcl']
    for executable in executables:
        copy(os.path.join(conf.BASEDIR, executable), jobdir)

    fp = open(os.path.join(jobdir, 'INPUT_INFO_STRUCTURES'), 'w')
    fp.write("%s %s\n" % (pdb1, pdb2))
    for chain in chains:
        resid1 = map(int, request.form.getlist('chain.%s.from' % chain))
        resid2 = map(int, request.form.getlist('chain.%s.to' % chain))
        fp.write("%s %s\n" % (chain, " ".join(['%d to %d' % (resid1[i], resid2[i]) for i in range(len(resid1))])))
    fp.close()

    fp = open(os.path.join(jobdir, 'run.com'), 'w')
    fc1 = request.form.get('fc1', 0.1)
    cutoff1 = request.form.get('cutoff1', 15.0)
    offset1 = request.form.get('offset1', 5.0)
    stepsize_cusp1 = request.form.get('stepsize_cusp1', 0.8)
    stepsize_slid1 = request.form.get('stepsize_slid1', 0.04)
    fc2 = request.form.get('fc2', 0.1)
    cutoff2 = request.form.get('cutoff2', 15.0)
    offset2 = request.form.get('offset2', 0.0)
    target = request.form.get('target_rmsd', 0.1)
    stepsize_cusp2 = request.form.get('stepsize_cusp2', 0.8)
    stepsize_slid2 = request.form.get('stepsize_slid2', 0.04)

    context = {}
    context['vmd'] = conf.VMD_EXECUTABLE
    for k in 'fc1 cutoff1 offset1 stepsize_cusp1 stepsize_slid1 fc2 cutoff2 offset2 target stepsize_cusp2 stepsize_slid2'.split():
        context[k] = locals()[k]
    fp.write(render_template('run.tpl', **context))

    fp = open(os.path.join(jobdir, 'run.pbs'), 'w')
    fp.write(render_template('run_pbs.tpl', uuid=jobid))
    fp.close()

    check_taskserver_status()
    submit_task(jobid, email)

    return redirect(url_for('success', uuid=jobid))

def allowed_file(filename):
    return '.' in filename and \
        filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

def _count_chain_resid(pdbfile):
    chains = {}
    for line in open(os.path.join(app.config['UPLOAD_FOLDER'], pdbfile)):
        if line.startswith("ATOM"):
            chainid = line[21:22]
            resid = int(line[22:26])
            name = line[12:17].strip()
            if name != 'CA': continue
            if not chains.has_key(chainid):
                chains[chainid] = [resid,resid]
            else:
                chains[chainid][1] = resid
    return chains

@app.route('/upload', methods=['GET', 'POST'])
def upload():
    if request.method == 'GET': return redirect(url_for('index'))

    error = {'code': 200}
    pdb1 = request.files['pdb1']
    pdb2 = request.files['pdb2']

    if not pdb1 and not pdb2:
        error={'code': 500, 'message': 'There was problem uploading PDB files'}
        return render_template('index.html', error=error)

    if bool(pdb1) != bool(pdb2):
        error={'code': 501, 'message': 'Please upload two PDB files'}
        return render_template('index.html', error=error)

    if not ( allowed_file(pdb1.filename) and allowed_file(pdb2.filename) ):
        error={'code': 502, 'message': 'Please upload PDB format files'}
        return render_template('index.html', error=error)

    filename = secure_filename(pdb1.filename)
    pdb1.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
    pdb1_chains = _count_chain_resid(filename)

    filename = secure_filename(pdb2.filename)
    pdb2.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
    pdb2_chains = _count_chain_resid(filename)

    if pdb1_chains != pdb2_chains:
        error={'code': 503, 'message': 'The chain lengths of two PDB are different'}
        return render_template('index.html', error=error)

    if filecmp.cmp(os.path.join(app.config['UPLOAD_FOLDER'], pdb1.filename), \
                   os.path.join(app.config['UPLOAD_FOLDER'], pdb2.filename)):
        error={'code': 504, 'message': 'Two PDB files are identical'}
        return render_template('index.html', error=error)        

    chains = []
    for k in sorted(pdb1_chains.keys()):
        chains.append((k, pdb1_chains[k]))

    return render_template('index.html', error=error, chains=chains, pdb1=pdb1, pdb2=pdb2, action=url_for('pathfinder'))

if __name__ == '__main__':
    app.debug = True
    app.run()
