import os
import pickle
from flask import Flask
from flask import (flash, current_app, session, redirect, render_template, request, url_for, send_from_directory)
from werkzeug.utils import secure_filename
from datetime import datetime
import logging
import traceback
import ast
import time

import run
from clipper import mail
from clipper.globals import *

app = Flask(__name__, instance_relative_config=True)
app.config.from_mapping(
    SECRET_KEY= "dev",
    UPLOAD_FOLDER= "upload",
    DOWNLOAD_FOLDER = result_folder_name,
    LOG_FOLDER = "log",
    DATA_FOLDER = "data",
    DATA_FILE = "data.pkl",
)
app.config.from_pyfile('config.py', silent=True)
app.add_url_rule('/', endpoint='index')

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
    return send_from_directory(uploads, filename, as_attachment=True)


@app.route('/<int:jobid>/result', methods=['GET'])
def download_results(jobid):
    arguments = create_arguments(jobid)
    print("\n\n\nin download_results()\n\n")

    try:
        run.main(arguments)
        output_name = arguments['output_name'] + '.zip'
        email = arguments['email']

        if email:
            try:
                mail.send_email(email, jobid, output_name)
            except Exception as err:
                print(f'Other error: {err}')
                traceback.print_exc()
                logging.info(f"Failed to send email for job id {session['jobid']}")

        downloads = arguments['resultdest']

        return send_from_directory(downloads, output_name, as_attachment=True)

    except TypeError as err:
        print(f'Typeerror: {err}')
        return redirect(url_for('error', jobid=jobid))

    except KeyError as err:
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
        session['form']['output_name'] = str(time.strftime("%Y%m%d%H%M%S")) + '_annotated_data'
    if 'email' not in session['form']:
        session['form']['email'] = False

   
    timestamp = datetime.now()
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
                    'output_filetype': session['form']['output_filetype'],
                    'separate': session['form']['separate'],
                    'pymol_verbose': False,
                    'jobid': str(jobid),
                    'logfile': logfile,
                    'resultdest': downloads,
                    'datafolder': datafolder,
                    'email': session['form']['email']
                    }
    except Exception as err:
        print(f'Other error: {err}')
        traceback.print_exc()

    return arguments
