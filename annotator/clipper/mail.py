import os
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

	current_app.config['MAIL_SERVER']='smtp.gmail.com'
	current_app.config['MAIL_PORT'] = 587
	current_app.config['MAIL_USERNAME'] = 'annotator.dtu@gmail.com'
	current_app.config['MAIL_PASSWORD'] = 'konka2612'
	current_app.config['MAIL_USE_TLS'] = True
	current_app.config['MAIL_USE_SSL'] = False

	mail = Mail(current_app)

	msg = Message(f'Annotator results {jobid}', sender = 'annotator.dtu@gmail.com', 
		recipients = [email_address])
	msg.body = "Thanks for using our tool. Your results are ready. \
You can find them as an atachment in this email."
	with current_app.open_resource(os.path.join(downloads, filename)) as fh:
		msg.attach(filename, filetype, fh.read())

	mail.send(msg)
	print(f'Sent result notification for {jobid}')
