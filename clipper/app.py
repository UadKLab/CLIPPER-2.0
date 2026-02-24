import logging
import os
import traceback

from flask import Flask
from flask import flash, redirect, render_template, request, send_from_directory, session, url_for

import run
from bin import mail
from bin.globals import result_folder_name
try:
    from .web import jobs
except ImportError:
    from web import jobs


app = Flask(__name__, instance_relative_config=True)
app.config.from_mapping(
    SECRET_KEY=os.getenv("CLIPPER_SECRET_KEY", "dev"),
    UPLOAD_FOLDER="upload",
    DOWNLOAD_FOLDER=result_folder_name,
    LOG_FOLDER="log",
    DATA_FOLDER="data",
    DATA_FILE="data.pkl",
)
app.config.from_pyfile("config.py", silent=True)
app.add_url_rule("/", endpoint="index")


@app.route("/", methods=["GET", "POST"])
def index():
    """Render the main submission page."""
    return render_template("index.html")


@app.route("/new_job", methods=["GET", "POST"])
def start_new_job():
    """Create a new job and store uploaded files/form values in session."""
    if request.method != "POST":
        return redirect(url_for("index"))

    jobs.ensure_runtime_files()

    try:
        jobid = jobs.reserve_job_id()
        jobs.store_submission(jobid)
    except ValueError as err:
        flash(str(err))
        return redirect(url_for("index"))
    except Exception as err:
        logging.exception("Failed to initialize new job")
        flash(f"Unable to start job: {err}")
        return redirect(url_for("index"))

    return redirect(url_for("submission", jobid=jobid))


@app.route("/<int:jobid>/submitted_job", methods=["GET"])
def submission(jobid):
    filename = session.get("infile")
    if not filename:
        flash("No active job found. Please resubmit your files.")
        return redirect(url_for("index"))

    return render_template("submission.html", jobid=jobid, filename=filename)


@app.route("/<int:jobid>/<filename>", methods=["GET"])
def download_input(jobid, filename):
    uploads = jobs.app_path(app.config["UPLOAD_FOLDER"])
    return send_from_directory(uploads, filename, as_attachment=True)


@app.route("/<int:jobid>/result", methods=["GET"])
def download_results(jobid):
    try:
        arguments = jobs.create_arguments(jobid)
        run.main(arguments)

        output_name = f"{arguments['output_name']}.zip"
        email = arguments["email"]

        if email:
            try:
                mail.send_email(email, jobid, output_name)
            except Exception:
                logging.exception("Failed to send email for job id %s", session.get("jobid"))

        return send_from_directory(arguments["resultdest"], output_name, as_attachment=True)

    except Exception as err:
        logging.exception("Job %s failed", jobid)
        session["error"] = str(err)
        session["traceback"] = traceback.format_exc()
        return redirect(url_for("error", jobid=jobid))


@app.route("/coming", methods=["GET"])
def coming():
    return render_template("coming_soon.html")


@app.route("/<int:jobid>/error", methods=["GET"])
def error(jobid):
    err = session.get("error", "Unknown error")
    return render_template("error.html", jobid=jobid, error=err)


if __name__ == "__main__":
    app.run(debug=True)
