from flask import Flask
from flask import render_template
app = Flask(__name__)

@app.errorhandler(404)
def page_not_found(e):
    return render_template('404.html'), 404

@app.route('/')
def main():
    return render_template('index.html')

@app.route('/pathfinder')
def pathfinder():
    return render_template('pathfinder.html')

if __name__ == '__main__':
    app.debug = True
    app.run()
