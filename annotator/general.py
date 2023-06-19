import os
import sys
import numpy as np
import pickle
from flask import (Blueprint, flash, current_app, session, redirect, render_template, request, url_for, send_from_directory)
from werkzeug.utils import secure_filename
from datetime import datetime
import logging
from pprint import pprint
import traceback
import ast

import os
from flask import Flask
from flask_mail import Mail

import run
from clipper import mail

#bp = Blueprint('general', __name__)

#def create_app(test_config=None):

#app = Flask(__name__, instance_relative_config=True, static_folder="clipper/statis", template_folder="clipper/template")
app = Flask(__name__, instance_relative_config=True)
app.config.from_mapping(
    SECRET_KEY= "dev",
    UPLOAD_FOLDER= "upload",
    DOWNLOAD_FOLDER = "download",
    LOG_FOLDER = "log",
    DATA_FOLDER = "data",
    DATA_FILE = "data.pkl",
)
app.config.from_pyfile('config.py', silent=True)
"""
if test_config is None:
    app.config.from_pyfile('config.py', silent=True)
else:
    app.config.from_mapping(test_config)
"""
#app.register_blueprint(bp)
app.add_url_rule('/', endpoint='index')

#return app

@app.route('/', methods=["GET", "POST"])
def index():
    if request.method == 'POST':

        data_folder = os.path.join(current_app.root_path, current_app.config['DATA_FOLDER'])
        data_path = os.path.join(data_folder, current_app.config['DATA_FILE'])
        upload_folder = os.path.join(current_app.root_path, current_app.config['UPLOAD_FOLDER'])

        with open(data_path, 'rb') as fh:
            dic = pickle.load(fh)
        
        jobid = dic['jobid']
        dic['jobid'] += 1
        
        with open(data_path, 'wb') as fh:
            pickle.dump(dic, fh)

        input_file = request.files['infile']
        cond_file = request.files['condfile']
        prot_file = request.files['protfile']

        if input_file and allowed_file(input_file.filename):
            input_filename = secure_filename(input_file.filename)
            input_file.save(os.path.join(upload_folder, input_filename))
            session['infile'] = input_filename

            if cond_file and allowed_file(cond_file.filename):
                cond_filename = secure_filename(cond_file.filename)
                cond_file.save(os.path.join(upload_folder, cond_filename))
                session['conditionfile'] = cond_filename
            else:
                session['conditionfile'] = None
            
            if prot_file and allowed_file(prot_file.filename):
                prot_filename = secure_filename(prot_file.filename)
                prot_file.save(os.path.join(upload_folder, prot_filename))
                session['proteasefile'] = prot_filename
            else:
                session['proteasefile'] = None

            session['jobid'] = jobid
            session['form'] = request.form

            return redirect(url_for('submission', jobid=jobid))

        else:
            flash('filename not valid')
            return(redirect(request.url))

    return render_template('index.html')


@app.route('/<int:jobid>/submit', methods=['GET'])
def submission(jobid):

    filename = session['infile']
    
    return render_template('submission.html', jobid=jobid, filename=filename)


@app.route('/<int:jobid>/<filename>', methods=['GET'])
def download_input(jobid, filename):
    uploads = os.path.join(current_app.root_path, current_app.config['UPLOAD_FOLDER'])
    
    return send_from_directory(directory=uploads, filename=filename)


@app.route('/<int:jobid>/result', methods=['GET'])
def download_results(jobid):
    arguments = create_arguments(jobid)
    print("in download_results()")
    try:
        outfile = run.main(arguments)
        outname = os.path.basename(outfile)
        print(f"outfile: {outfile}")
        print(f'putname: {outname}')
        email = arguments['email']
        print(f'email: {email}')

        if email:
            try:
                mail.send_email(email, jobid, outname)
            except:
                logging.info(f"Failed to send email for {session['jobid']}")

        downloads = arguments['resultdest']
        print('downloads')
        print(downloads)
        return send_from_directory(directory=downloads, filename=outname, as_attachment=True)

    except TypeError as err:
        print(f'Typeerror: {err}')
        session['error'] = str(err)
        return redirect(url_for('error', jobid=jobid))

    except KeyError as err:
        session['error'] = str(sys.exc_info()[0])
        print(f"Key error: {err}")
        traceback.print_exc()
        return redirect(url_for('error', jobid=jobid))
    except Exception as err:
        print(f'Other error: {err}')
        traceback.print_exc()
        return redirect(url_for('error', jobid=jobid))


@app.route('/coming', methods=['GET'])
def coming():

    return render_template('coming_soon.html')


@app.route('/<int:jobid>/error', methods=['GET'])
def error(jobid):

    print('\nIn error()\n')

    return render_template('error.html', jobid=jobid, error=session['error'])


def allowed_file(filename):
    ALLOWED_EXTENSIONS = {'txt', 'xlsx', 'xls', 'csv'}
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


def create_arguments(jobid):

    uploads = os.path.join(current_app.root_path, current_app.config['UPLOAD_FOLDER'])
    downloads = os.path.join(current_app.root_path, current_app.config['DOWNLOAD_FOLDER'])
    logfile = os.path.join(current_app.root_path, current_app.config['LOG_FOLDER'], str(jobid) + '.log')
    datafolder = os.path.join(current_app.root_path, current_app.config['DATA_FOLDER'])
    email = False
    cond_file = None
    prot_file = None

    if session['conditionfile']:
        cond_file = os.path.join(uploads, session['conditionfile'])
    if session['proteasefile']:
        prot_file = os.path.join(uploads, session['proteasefile'])
    if 'statistic' not in session['form']:
        session['form']['statistic'] = False
    if 'stat_pairwise' not in session['form']:
        session['form']['stat_pairwise'] = False
    if 'significance' not in session['form']:
        session['form']['significance'] = None
    if 'dropna' not in session['form']:
        session['form']['dropna'] = False
    if 'fillna' not in session['form']:
        session['form']['fillna'] = None
    if 'singlecpu' not in session['form']:
        session['form']['singlecpu'] = False
    if 'noexopeptidase' not in session['form']:
        session['form']['noexopeptidase'] = False
    if 'nomerops' not in session['form']:
        session['form']['nomerops'] = False
    if 'visualize' not in session['form']:
        session['form']['visualize'] = False
    if 'enrichment' not in session['form']:
        session['form']['enrichment'] = False
    if 'pathway' not in session['form']:
        session['form']['pathway'] = False
    if 'logo' not in session['form']:
        session['form']['logo'] = None
    if 'pseudocounts' not in session['form']:
        session['form']['pseudocounts'] = False
    if 'separate' not in session['form']:
        session['form']['separate'] = False
    if 'output_name' not in session['form']:
        session['form']['output_name'] = None
    if 'email' not in session['form']:
        email = session['form']['email'] = None

   
    timestamp = datetime.now()
    ftimestamp = timestamp.strftime(format="%d%m%y%H%M%S")
    formatted_timestamp = timestamp.strftime(format="%A %B %d %Y, %H:%M:%S")
    logging.basicConfig(filename=logfile, filemode="w", level=logging.INFO)
    logging.info(f"Annotator started, {formatted_timestamp}")
  
    
    try:
        arguments = {
                    'infile': os.path.join(uploads, session['infile']),
                    'infile_type': session['form']['infile_type'],
                    'software': session['form']['software'],
                    'level': session['form']['filter'],
                    'dropna': session['form']['dropna'],
                    'fillna': session['form']['fillna'],
                    'sleeptime': float(session['form']['sleeptime']),
                    'noexo': session['form']['noexopeptidase'],
                    'nomerops': session['form']['nomerops'],
                    'calcstructure': ast.literal_eval(session['form']['calcstructure']),
                    'singlecpu': session['form']['singlecpu'],
                    'conditionfile': cond_file,
                    'proteasefile': prot_file,
                    'stat': session['form']['statistic'], 
                    'stat_pairwise': session['form']['stat_pairwise'], 
                    'significance': session['form']['significance'],
                    'visualize': session['form']['visualize'],
                    'logo': session['form']['logo'],
                    'pseudocounts': session['form']['pseudocounts'],
                    'cleavagevis': session['form']['cleavagevis'],
                    'enrichment': session['form']['enrichment'],
                    'pathway': session['form']['pathway'],
                    'output_name': session['form']['output_name'],
                    'outfile_type': session['form']['outfile_type'],
                    'separate': session['form']['separate'],
                    'pymol_verbose': False,
                    'jobid': str(jobid),
                    'logfile': logfile,
                    'resultdest': downloads,
                    'datafolder': datafolder,
                    'email': email
                    }
    except Exception as err:
        print(f'ERROR: {err}')


    return arguments
