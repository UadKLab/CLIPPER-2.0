from __future__ import annotations

import logging
import pickle
import time
from datetime import datetime
from pathlib import Path

from flask import current_app, request, session
from werkzeug.utils import secure_filename


DEFAULT_OUTPUT_SUFFIX = "_annotated_data"
ALLOWED_EXTENSIONS = {"txt", "xlsx", "xls", "csv"}
TRUTHY_VALUES = {"1", "true", "yes", "on"}


def app_path(*parts: str) -> Path:
    return Path(current_app.root_path, *parts)


def allowed_file(filename: str) -> bool:
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


def coerce_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in TRUTHY_VALUES


def read_form_bool(form_data: dict, *keys: str) -> bool:
    for key in keys:
        if key in form_data:
            return coerce_bool(form_data.get(key))
    return False


def normalize_optional(value):
    if value is None:
        return None
    value = str(value).strip()
    return value or None


def ensure_runtime_files() -> None:
    upload_dir = app_path(current_app.config["UPLOAD_FOLDER"])
    data_dir = app_path(current_app.config["DATA_FOLDER"])
    log_dir = app_path(current_app.config["LOG_FOLDER"])
    result_dir = app_path(current_app.config["DOWNLOAD_FOLDER"])
    data_file = data_dir / current_app.config["DATA_FILE"]

    for folder in (upload_dir, data_dir, log_dir, result_dir):
        folder.mkdir(parents=True, exist_ok=True)

    if not data_file.exists():
        with data_file.open("wb") as fh:
            pickle.dump({"jobid": 1}, fh)


def reserve_job_id() -> int:
    data_path = app_path(current_app.config["DATA_FOLDER"], current_app.config["DATA_FILE"])
    with data_path.open("rb") as fh:
        payload = pickle.load(fh)

    jobid = int(payload.get("jobid", 1))
    payload["jobid"] = jobid + 1

    with data_path.open("wb") as fh:
        pickle.dump(payload, fh)

    return jobid


def save_optional_upload(file_obj, upload_folder: Path):
    if not file_obj or not file_obj.filename:
        return None
    if not allowed_file(file_obj.filename):
        return None

    filename = secure_filename(file_obj.filename)
    if not filename:
        return None

    destination = upload_folder / filename
    file_obj.save(destination)
    return filename


def store_submission(jobid: int) -> None:
    upload_folder = app_path(current_app.config["UPLOAD_FOLDER"])

    input_file = request.files.get("infile")
    if input_file is None or not input_file.filename or not allowed_file(input_file.filename):
        raise ValueError("Input file is missing or has an unsupported file extension")

    input_filename = secure_filename(input_file.filename)
    input_file.save(upload_folder / input_filename)

    session["jobid"] = jobid
    session["infile"] = input_filename
    session["conditionfile"] = save_optional_upload(request.files.get("condfile"), upload_folder)
    session["proteasefile"] = save_optional_upload(request.files.get("protfile"), upload_folder)
    session["form"] = request.form.to_dict(flat=True)


def default_output_name() -> str:
    return f"{time.strftime('%Y%m%d%H%M%S')}{DEFAULT_OUTPUT_SUFFIX}"


def create_arguments(jobid):
    uploads = app_path(current_app.config["UPLOAD_FOLDER"])
    downloads = app_path(current_app.config["DOWNLOAD_FOLDER"])
    logfile = app_path(current_app.config["LOG_FOLDER"], f"{jobid}.log")
    datafolder = app_path(current_app.config["DATA_FOLDER"])

    form = dict(session.get("form") or {})

    infile = session.get("infile")
    if not infile:
        raise ValueError("Input file is missing from session")

    cond_file = None
    prot_file = None
    if session.get("conditionfile"):
        cond_file = str(uploads / session["conditionfile"])
    if session.get("proteasefile"):
        prot_file = str(uploads / session["proteasefile"])

    calcstructure = normalize_optional(form.get("calcstructure"))
    cleavagevis = normalize_optional(form.get("cleavagevis"))
    logo = normalize_optional(form.get("logo"))
    significance = normalize_optional(form.get("significance"))
    output_name = normalize_optional(form.get("outfile_name")) or default_output_name()
    email = normalize_optional(form.get("email"))

    logging.basicConfig(filename=logfile, filemode="w", level=logging.INFO)
    formatted_timestamp = datetime.now().strftime("%A %B %d %Y, %H:%M:%S")
    logging.info("Annotator started, %s", formatted_timestamp)

    arguments = {
        "infile": str(uploads / infile),
        "infile_type": form.get("infile_type", "infer"),
        "software": form.get("software", "infer"),
        "level": form.get("filter", "all"),
        "dropna": read_form_bool(form, "dropna"),
        "fillna": normalize_optional(form.get("fillna")),
        "alpha": float(form.get("alpha", 0.05)),
        "sleeptime": float(form.get("sleeptime", 0.2)),
        "noexo": read_form_bool(form, "noexopeptidase"),
        "nomerops": read_form_bool(form, "nomerops"),
        "calcstructure": calcstructure,
        "threadingcores": "max",
        "conditionfile": cond_file,
        "proteasefile": prot_file,
        "stat": read_form_bool(form, "statistics", "statistic"),
        "stat_pairwise": read_form_bool(form, "stat_pairwise"),
        "significance": significance,
        "multipletesting": read_form_bool(form, "multipletesting"),
        "multipletestingmethod": "fdr_bh",
        "visualize": read_form_bool(form, "visualize"),
        "logo": logo,
        "logo_fc": float(form.get("logo_fc", 3)),
        "volcano_foldchange": float(form.get("volcano_foldchange", 1.5)),
        "cleavagesitesize": int(form.get("cleavagesitesize", 4)),
        "pseudocounts": read_form_bool(form, "pseudocounts"),
        "cleavagevis": cleavagevis,
        "enrichment": read_form_bool(form, "enrichment"),
        "pathway": read_form_bool(form, "pathway"),
        "output_name": output_name,
        "output_filetype": form.get("output_filetype", "xlsx"),
        "separate": read_form_bool(form, "separate"),
        "pymol_verbose": False,
        "jobid": str(jobid),
        "logfile": str(logfile),
        "resultdest": str(downloads),
        "datafolder": str(datafolder),
        "email": email,
    }

    return arguments
