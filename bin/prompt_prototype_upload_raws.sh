#!/bin/bash
# This file is part of prompt_prototype.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# This script uploads the raw files from the HSC PDR2 run to Google Storage. It
# renames the files to match prompt_prototype conventions. The user must have
# gsutil already configured.

set -e  # Abort on any error

RAW_DIR="/datasets/hsc/raw/ssp_pdr2/2016-03-07"
UPLOAD_BUCKET=rubin-prompt-proto-unobserved

# Filename format is defined in activator/raw.py:
# instrument/detector/group/snap/exposureId/filter/instrument-group-snap-exposureId-filter-detector
gsutil cp "${RAW_DIR}/HSCA05913553.fits" \
    gs://${UPLOAD_BUCKET}/HSC/0/2016030700001/0/0059134/HSC-G/HSC-2016030700001-0-0059134-HSC-G-0.fits
gsutil cp "${RAW_DIR}/HSCA05913542.fits" \
    gs://${UPLOAD_BUCKET}/HSC/4/2016030700001/0/0059134/HSC-G/HSC-2016030700001-0-0059134-HSC-G-4.fits
gsutil cp "${RAW_DIR}/HSCA05913543.fits" \
    gs://${UPLOAD_BUCKET}/HSC/5/2016030700001/0/0059134/HSC-G/HSC-2016030700001-0-0059134-HSC-G-5.fits

gsutil cp "${RAW_DIR}/HSCA05914353.fits" \
    gs://${UPLOAD_BUCKET}/HSC/0/2016030700002/0/0059142/HSC-G/HSC-2016030700002-0-0059142-HSC-G-0.fits
gsutil cp "${RAW_DIR}/HSCA05914343.fits" \
    gs://${UPLOAD_BUCKET}/HSC/5/2016030700002/0/0059142/HSC-G/HSC-2016030700002-0-0059142-HSC-G-5.fits
gsutil cp "${RAW_DIR}/HSCA05914337.fits" \
    gs://${UPLOAD_BUCKET}/HSC/11/2016030700002/0/0059142/HSC-G/HSC-2016030700002-0-0059142-HSC-G-11.fits

gsutil cp "${RAW_DIR}/HSCA05915112.fits" \
    gs://${UPLOAD_BUCKET}/HSC/50/2016030700003/0/0059150/HSC-G/HSC-2016030700003-0-0059150-HSC-G-50.fits
gsutil cp "${RAW_DIR}/HSCA05915116.fits" \
    gs://${UPLOAD_BUCKET}/HSC/58/2016030700003/0/0059150/HSC-G/HSC-2016030700003-0-0059150-HSC-G-58.fits

gsutil cp "${RAW_DIR}/HSCA05916109.fits" \
    gs://${UPLOAD_BUCKET}/HSC/43/2016030700004/0/0059160/HSC-G/HSC-2016030700004-0-0059160-HSC-G-43.fits
gsutil cp "${RAW_DIR}/HSCA05916113.fits" \
    gs://${UPLOAD_BUCKET}/HSC/51/2016030700004/0/0059160/HSC-G/HSC-2016030700004-0-0059160-HSC-G-51.fits
