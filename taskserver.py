import socket
import threading
import SocketServer

from flask import g, Flask, render_template
import os
from main import get_db, get_job_folder, download
import conf
import errno
from flask_mail import Mail, Message

NUMJOBS_CONCUR = 4

def pid_exists(pid):
    try:
        os.kill(pid, 0)
    except OSError, e:
        return e.errno == errno.EPERM
    else:
        return True

def check_queue(ip, port):
    app = Flask(__name__)
    app.config['SERVER_NAME'] = conf.SERVER_NAME
    app.add_url_rule('/pathfinder/download/<uuid>', 'download', download)

    with app.app_context():
        db = get_db()
        cur = db.cursor()
        cur.execute('select * from job where status == ?', ('R',))
        numjobs = 0
        mail = None
        for row in cur.fetchall():
            uuid, email, date, status = row
            # finished?
            pid = int(open(os.path.join(get_job_folder(uuid), 'run.pid'), 'r').read())
            if pid_exists(pid):
                numjobs += 1
            else:
                mail = Mail(app)
                msg = Message('[PathFinder] Your job is finished.', sender='roux@uchicago.edu', recipients=[email, 'sunhwanj@gmail.com'])

                if os.path.exists(os.path.join(get_job_folder(uuid), 'pathway.pdb')):
                    cur.execute('update job set status = ? where uuid = ?', ('F',uuid))
                    msg.body = render_template('email.tpl', uuid=uuid, has_error=False)
                    msg.attach('pathway.pdb', 'text/pdb', open(os.path.join(get_job_folder(uuid), 'pathway.pdb')).read())
                else:
                    cur.execute('update job set status = ? where uuid = ?', ('E',uuid))
                    msg.body = render_template('email.tpl', uuid=uuid, has_error=True)

        if numjobs < NUMJOBS_CONCUR:
            cur.execute('select * from job where status == ?', ('Q',))
            for row in cur.fetchall():
                uuid, email, date, status = row
                newpid = client(ip, port, "SPAWN:%s" % uuid)
                open(os.path.join(get_job_folder(uuid), 'run.pid'), 'w').write(newpid)
                cur.execute('update job set status = ? where uuid = ?', ('R',uuid))
                numjobs += 1
                if numjobs >= NUMJOBS_CONCUR: break

        db.commit()
        db.close()

        if mail: mail.send(msg)


class ThreadedTCPRequestHandler(SocketServer.BaseRequestHandler):

    def handle(self):
        data = self.request.recv(1024)
        handler, data = data.split(":")
        #cur_thread = threading.current_thread()
        func = getattr(self, 'handle_%s' % handler.lower())
        newpid = os.fork()
        if newpid == 0:
            func(data)
        else:
            self.request.sendall(str(newpid))

    def handle_spawn(self, uuid):
        cwd = os.getcwd()
        os.chdir(get_job_folder(uuid))
        os.system('/bin/csh ./run.com > run.out')


class ThreadedTCPServer(SocketServer.ForkingMixIn, SocketServer.TCPServer):
    pass


def client(ip, port, message):
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect((ip, port))
    try:
        sock.sendall(message)
        response = sock.recv(1024)
        return response
    finally:
        sock.close()

if __name__ == "__main__":
    # Port 0 means to select an arbitrary unused port
    HOST, PORT = "localhost", 0

    server = ThreadedTCPServer((HOST, PORT), ThreadedTCPRequestHandler)
    ip, port = server.server_address
    open('taskserver.info', 'w').write("%s:%s" % (ip, port))

    # Start a thread with the server -- that thread will then start one
    # more thread for each request
    server_thread = threading.Thread(target=server.serve_forever)
    # Exit the server thread when the main thread terminates
    server_thread.daemon = True
    server_thread.start()

    import time
    while 1:
        check_queue(ip, port)
        time.sleep(10)

    print "Server loop running in thread:", server_thread.name
