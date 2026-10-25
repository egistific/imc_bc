import os
import shutil
import yaml
import glob
import pandas as pd


configfile: "config.yaml"

ROOT_DIR = config["project"]["root_dir"]
SEED = config["project"]["seed"]
RAW_DIR = config["data"]["raw_dir"]
METADATA_DIR = config["data"]["metadata_dir"]
STEINBOCK_DIR = config["data"]["steinbock_dir"]

RAW_FILES = glob.glob(f"{RAW_DIR}/*.mcd", recursive=True) + glob.glob(f"{RAW_DIR}/*.txt", recursive=True) + glob.glob(f"{RAW_DIR}/*.zip", recursive=True)
PATIENT_FILE = os.path.join(METADATA_DIR, config["data"]["patient_file"])
SM_FILE = os.path.join(METADATA_DIR, config["data"]["compensation_file"])
MARKERS_YAML = os.path.join(METADATA_DIR, config["analysis"]["astir"]["markers_yaml"])

PANEL = os.path.join(STEINBOCK_DIR, config["panel"]["original"])
PANEL_EDITED = os.path.join(STEINBOCK_DIR, config["panel"]["edited"])
IMG_METADATA = os.path.join(STEINBOCK_DIR, config["preprocessing"]["img_metadata"])

CHANNELS_DISCARD = config["channels"]["discard"]
CHANNELS_NUCLEAR = config["channels"]["nuclear"]
CHANNELS_MEMBRANE_OR_CYTOPLASM = config["channels"]["membrane_or_cytoplasm"]
HPF = config["preprocessing"]["hot_pixel_filter"]
IMG_DIR = config["preprocessing"]["img_dir"]
METHOD = config["cell_segmentation"]["method"]
APP = config["cell_segmentation"]["app"]
IS_MINMAX = bool(config["cell_segmentation"]["minmax"])
MASKS_DIR = config["cell_segmentation"]["masks_dir"]
IMG_WITH_MASKS_DIR = config["img_with_masks_dir"]
INTENSITIES_DIR = config["measure_intensities"]["intensities_dir"]
REGIONPROPS_DIR = config["measure_regionprops"]["regionprops_dir"]
NEIGHBORS_DIR = config["measure_neighbors"]["neighbors_dir"]
EXPANSION = config["measure_neighbors"]["expansion"]
DMAX = config["measure_neighbors"]["dmax"]
DNA1_CUTOFF = config["analysis"]["dna1_cutoff"]
ANALYSIS_DIR = config["analysis"]["output_dir"]
RDS_FILE = os.path.join(ANALYSIS_DIR, config["analysis"]["rds_file"])
RDS_FILE_COMPENSATED = os.path.join(ANALYSIS_DIR, config["analysis"]["rds_file_compensated"])
RDS_FILE_FILTERED = os.path.join(ANALYSIS_DIR, config["analysis"]["rds_file_filtered"])
EXPRS_FILE = os.path.join(ANALYSIS_DIR, config["analysis"]["astir"]["exprs_file"])
DESIGN_FILE = os.path.join(ANALYSIS_DIR, config["analysis"]["astir"]["design_file"])
BATCH_CORRECTION = config["analysis"]["astir"]["batch_correction"]
CELLTYPES_FILE = os.path.join(ANALYSIS_DIR, config["analysis"]["astir"]["celltypes_file"])
FIGURES_DIR = os.path.join(ANALYSIS_DIR, config["analysis"]["plots"]["figures_dir"])
MAX_EPOCHS = config["analysis"]["astir"]["max_epochs"]
BATCH_SIZE = config["analysis"]["astir"]["batch_size"]
LEARNING_RATE = config["analysis"]["astir"]["learning_rate"]
N_INIT = config["analysis"]["astir"]["num_init"]
N_INIT_EPOCHS = config["analysis"]["astir"]["num_init_epochs"]
DELTA_LOSS = config["analysis"]["astir"]["delta_loss"]
THRESHOLD = config["analysis"]["astir"]["threshold"]


def run_docker(root_dir, args):
    return f"docker run --rm -v {root_dir}:/data -w /data -u $(id -u):$(id -g) ghcr.io/bodenmillergroup/steinbock:0.16.1 {args}"

def panel_args(raw_dir, panel):
    return f"preprocess imc panel --mcd {raw_dir} --txt {raw_dir} --unzip -o {panel}"

def new_panel(old_filename, new_filename, method, channels_discard, channels_nuclear, channels_membrane_or_cytoplasm):
    df = pd.read_csv(old_filename)
    print(df)
    if method not in df.columns:
        df[method] = pd.NA
    df.loc[df["channel"].isin(channels_discard), "keep"] = 0
    df.loc[df["channel"].isin(channels_nuclear), method] = 1
    df.loc[df["channel"].isin(channels_membrane_or_cytoplasm), method] = 2
    df[method] = pd.to_numeric(df[method], errors="coerce").replace([float("inf"), float("-inf")], pd.NA).astype("Int64")
    df.to_csv(new_filename, index=False)

def preprocess_args(panel, hpf, img_dir, img_metadata):
    return f"preprocess imc images --panel {panel} --hpf {hpf} --imgout {img_dir} --infoout {img_metadata}"

def segment_args(method, app, img_dir, panel, is_minmax, out_dir):
    cmd = f"segment {method} --app {app} --panel {panel} --img {img_dir}"
    if is_minmax:
        cmd += " --minmax"
    cmd += f" -o {out_dir}"
    print(cmd)
    return cmd

def images_with_masks(img_dir, masks_dir, img_with_masks_dir):
    image_paths =glob.glob(os.path.join(img_dir, "*.tiff"))
    keep_paths = [img for img in image_paths if os.path.isfile(os.path.join(masks_dir, os.path.basename(img)))]
    if not os.path.exists(img_with_masks_dir):
        os.makedirs(img_with_masks_dir)
    for source_path in keep_paths:
        dest_path = os.path.join(img_with_masks_dir, os.path.basename(source_path))
        print(dest_path)
        shutil.copy2(source_path, dest_path)

def intensities_arg(img_dir, masks_dir, panel):
    return f"measure intensities --img {img_dir} --masks {masks_dir} --panel {panel}"

def regionprops_args(img_dir, masks_dir):
    return f"measure regionprops --img {img_dir} --masks {masks_dir}"

def neighbors_args(masks_dir, expansion, dmax, neighbors_dir):
    return f"measure neighbors --masks {masks_dir} --type {expansion} --dmax {dmax} -o {neighbors_dir}"


rule all:
    input:
        PANEL,
        PANEL_EDITED,
        IMG_DIR,
        IMG_METADATA,
        MASKS_DIR,
        INTENSITIES_DIR,
        REGIONPROPS_DIR,
        NEIGHBORS_DIR,
        RDS_FILE,
        RDS_FILE_COMPENSATED,
        RDS_FILE_FILTERED,
        CELLTYPES_FILE
        # COMBINED_DIR,
        # COMBINED_FILENAME

rule create_panel:
    input:
        raw=RAW_FILES
    output:
        panel=PANEL
    params:
        cmd=lambda wildcards: run_docker(ROOT_DIR, panel_args(RAW_DIR, PANEL))
    log:
        "logs/create_panel.log"
    shell:
        """
        {params.cmd} > {log} 2>&1
        echo "Edit panel.csv and save as panel_edited.csv."
        """

rule edit_panel:
    input:
        panel=PANEL
    output:
        panel=PANEL_EDITED
    log:
        "logs/edit_panel.log"
    run:
        new_panel(input.panel, output.panel, METHOD, CHANNELS_DISCARD, CHANNELS_NUCLEAR, CHANNELS_MEMBRANE_OR_CYTOPLASM)

rule preprocess_images:
    input:
        panel=PANEL_EDITED,
        raw=RAW_FILES
    output:
        images=directory(IMG_DIR),
        info=IMG_METADATA
    params:
        cmd=lambda wildcards: run_docker(ROOT_DIR, preprocess_args(PANEL_EDITED, HPF, IMG_DIR, IMG_METADATA))
    log:
        "logs/preprocess_images.log"
    shell:
        """
         {params.cmd} > {log} 2>&1
        """

rule segment_cells:
    input:
        images=IMG_DIR, # directory not valid for inputs
        panel=PANEL_EDITED,
        # info=IMG_METADATA
    output:
        masks=directory(MASKS_DIR)
    params:
        cmd=lambda wildcards: run_docker(ROOT_DIR, segment_args(METHOD, APP, IMG_DIR, PANEL_EDITED, IS_MINMAX, MASKS_DIR))
    log:
        "logs/segment_cells.log"
    shell:
        """
        {params.cmd} > {log} 2>&1
        """

rule images_with_masks:
    input:
        images=IMG_DIR,
        masks=MASKS_DIR
    output:
        images=directory(IMG_WITH_MASKS_DIR)
    log:
        "logs/img_with_masks.log"
    run:
        images_with_masks(input.images, input.masks, output.images)

rule measure_intensities:
    input:
        images=IMG_WITH_MASKS_DIR, # directory not valid for inputs
        masks=MASKS_DIR,
        panel=PANEL_EDITED,
    output:
        intensities=directory(INTENSITIES_DIR)
    params:
        cmd=lambda wildcards: run_docker(ROOT_DIR, intensities_arg(IMG_WITH_MASKS_DIR, MASKS_DIR, PANEL_EDITED))
    log:
        "logs/measure_intensities.log"
    shell:
        """
        {params.cmd} > {log} 2>&1
        """

rule measure_regionprops:
    input:
        images=IMG_WITH_MASKS_DIR, # directory not valid for inputs
        masks=MASKS_DIR,
        # panel=PANEL_EDITED,
    output:
        regionprops=directory(REGIONPROPS_DIR)
    params:
        cmd=lambda wildcards: run_docker(ROOT_DIR, regionprops_args(IMG_WITH_MASKS_DIR, MASKS_DIR))
    log:
        "logs/measure_regionprops.log"
    shell:
        """
        {params.cmd} > {log} 2>&1
        """

rule measure_neighbors:
    input:
        masks=MASKS_DIR
    output:
        neighbors=directory(NEIGHBORS_DIR)
    params:
        cmd=lambda wildcards: run_docker(ROOT_DIR, neighbors_args(MASKS_DIR, EXPANSION, DMAX, NEIGHBORS_DIR))
    log:
        "logs/measure_neighbors.log"
    shell:
        """
        {params.cmd} > {log} 2>&1
        """

rule import_steinbock:
    input:
        intensities=INTENSITIES_DIR,
        regionprops=REGIONPROPS_DIR,
        neighbors=NEIGHBORS_DIR,
        panel=PANEL_EDITED,
        patient=PATIENT_FILE,
        img_info=IMG_METADATA
    output:
        spe=RDS_FILE
    log:
        "logs/import_steinbock.log"
    script:
        "scripts/01_import_steinbock.R"

rule compensate_spillover:
    input:
        spe=RDS_FILE,
        sm=SM_FILE
    output:
        spe=RDS_FILE_COMPENSATED
    log:
        "logs/compensate_spillover.log"
    script:
        "scripts/02_compensate_spillover.R"

rule filter_cells:
    input:
        spe=RDS_FILE_COMPENSATED,
    output:
        spe=RDS_FILE_FILTERED,
        exprs=EXPRS_FILE,
        design=DESIGN_FILE
    log:
        "logs/filter_cells.log"
    script:
        "scripts/03_filter_cells.R"

rule run_astir:
    input:
        exprs=EXPRS_FILE,
        design=DESIGN_FILE,
        yaml=MARKERS_YAML,
    output:
        celltypes=CELLTYPES_FILE
    params:
        design_arg = f" --design {DESIGN_FILE}" if BATCH_CORRECTION is not None else "",
        seed=SEED,
        max_epochs=MAX_EPOCHS,
        batch_size=BATCH_SIZE,
        learning_rate=LEARNING_RATE,
        n_init=N_INIT,
        n_init_epochs=N_INIT_EPOCHS,
        delta_loss=DELTA_LOSS,
        threshold=THRESHOLD
    log:
        "logs/run_astir.log"
    shell:
        "python scripts/04_run_astir.py --exprs {input.exprs}{params.design_arg} --markers {input.yaml} --output {output.celltypes} \
        --seed {params.seed} --max_epochs {params.max_epochs} --batch_size {params.batch_size} --learning_rate {params.learning_rate} \
        --n_init {params.n_init} --n_init_epochs {params.n_init_epochs} --delta_loss {params.delta_loss} --threshold {params.threshold} > {log} 2>&1"

"""
rule plot_celltypes:
    input:
        spe=RDS_FILE_FILTERED,
        celltypes=CELLTYPES_FILE
    output:
        heatmap_mean=os.path.join(FIGURES_DIR, "heatmap_celltype_mean.pdf"),
        heatmap_cells=os.path.join(FIGURES_DIR, "heatmap_single_cells.pdf"),
        max_prob=os.path.join(FIGURES_DIR, "astir_max_probability.pdf"),
        composition=os.path.join(FIGURES_DIR, "composition_by_tissue.pdf"),
        image_props=os.path.join(FIGURES_DIR, "proportions_per_image.pdf"),
        spatial=os.path.join(FIGURES_DIR, "spatial_example_images.pdf"),
        props_csv=os.path.join(ROOT_DIR, ANALYSIS_DIR, "celltype_proportions_per_image.csv")
    log:
        "logs/plot_celltypes.log"
    script:
        "scripts/05_plot_celltypes.R"
"""