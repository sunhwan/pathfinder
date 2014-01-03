import sys
sys.path.append('/home/sunhwan/.local/lib/python2.7/site-packages/')
sys.path.append('/home/sunhwan/www/anmpathway/pathfinder/')

from flask import g, Flask, render_template
import os
from main import get_db, get_job_folder, download, status
import conf
import errno
from flask_mail import Mail, Message
import time

uuid = sys.argv[1]

if __name__ == "__main__":
    app = Flask(__name__, template_folder=os.path.join(conf.BASEDIR, 'templates'))
    app.config['SERVER_NAME'] = conf.SERVER_NAME
    app.add_url_rule('/anmpathway.cgi/download/<uuid>/<filename>', 'download', download)
    app.add_url_rule('/anmpathway.cgi/download/<uuid>', 'download', download)
    app.add_url_rule('/anmpathway.cgi/status/<uuid>', 'status', status)

    mail = None

    with app.app_context():
        db = get_db()
        cur = db.cursor()
        cur.execute('select * from job where uuid == ?', (uuid,))
        uuid, email, date, status = cur.fetchone()

        open(os.path.join(get_job_folder(uuid), 'run.pid'), 'w').write(str(os.getpid()))
        cur.execute('update job set status = ? where uuid = ?', ('R',uuid))
        db.commit()

        mail = Mail(app)
        msg = Message('[ANMPathway] Your job has started.', sender=conf.ADMIN_EMAIL, recipients=[email], bcc=[conf.ADMIN_EMAIL, conf.EXTRA_EMAIL])
        msg.html = render_template('email_started.tpl', uuid=uuid)
        mail.send(msg)

        etime = time.time()
        cwd = os.getcwd()
        os.chdir(get_job_folder(uuid))
        os.system('/bin/csh ./run.com > run.out')
        etime = time.time() - etime

        mail = Mail(app)
        msg = Message('[ANMPathway] Your job is finished.', sender=conf.ADMIN_EMAIL, recipients=[email], bcc=[conf.ADMIN_EMAIL, conf.EXTRA_EMAIL])


        if os.path.exists(os.path.join(get_job_folder(uuid), 'pathway.pdb')):
            cur.execute('update job set status = ? where uuid = ?', ('F',uuid))
            msg.html = render_template('email.tpl', uuid=uuid, has_error=False, etime=etime)
            #msg.attach('pathway.pdb', 'text/pdb', open(os.path.join(get_job_folder(uuid), 'pathway.pdb')).read())
            #msg.attach('close_contacts_5.0_10.0.txt', 'text/plain', open(os.path.join(get_job_folder(uuid), 'close_contacts_5.0_10.0')).read())

            close_contact = os.path.join(get_job_folder(uuid), 'close_contacts_5.0_10.0')
            if os.path.exists(close_contact) and os.stat(close_contact).st_size < 1:
                open(close_contact, 'w').write('no close residue pairs were found')
        else:
            cur.execute('update job set status = ? where uuid = ?', ('E',uuid))
            msg.html = render_template('email.tpl', uuid=uuid, has_error=True, etime=etime)

        db.commit()
        db.close()

        if mail:
            mail.send(msg)

