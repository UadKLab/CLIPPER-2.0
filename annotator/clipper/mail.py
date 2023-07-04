import os
import json
from flask import current_app, send_from_directory
from flask_mail import Message, Mail

mimetype = {'json' : 'application/json',
			'xlsx' : 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
			'xls' : 'application/vnd.ms-excel',
			'csv' : 'text/csv',
			'tsv' : 'text/tab-separated-values',
			'zip' : 'application/zip',
			}

def send_email(email_address, jobid, filename):

	downloads = os.path.join(current_app.root_path, current_app.config['DOWNLOAD_FOLDER'])
	extension = filename.rsplit('.', 1)[1]
	filetype = mimetype[extension]

	print(os.getcwd())
	with open('data/credentials.json', 'r') as file:
		credentials = json.load(file)

	email = credentials['email']
	password = credentials['password']

	current_app.config['MAIL_SERVER']='smtp.gmail.com'
	current_app.config['MAIL_PORT'] = 587
	current_app.config['MAIL_USERNAME'] = email
	current_app.config['MAIL_PASSWORD'] = password
	current_app.config['MAIL_USE_TLS'] = True
	current_app.config['MAIL_USE_SSL'] = False

	mail = Mail(current_app)

	msg = Message(f'Annotator results {jobid}', sender = email, 
		recipients = [email_address])
	msg.body = "Thanks for using our tool. Your results are ready. You can find them as an atachment in this email."
	with current_app.open_resource(os.path.join(downloads, filename)) as fh:
		msg.attach(filename, filetype, fh.read())

	mail.send(msg)
	print(f'Sent result notification for {jobid}')
