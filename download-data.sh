
set -e
set -o pipefail

# Use the OpenPBTA bucket as the default.
URL=${OPENPEDCAN_URL:-https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/pbta-alt-2022}
RELEASE=${OPENPEDCAN_RELEASE:-v1}
PREVIOUS=${OPENPEDCAN_RELEASE:-NA}

# Remove old symlinks in data
find data -type l -delete

# The md5sum file provides our single point of truth for which files are in a release.
curl --create-dirs $URL/$RELEASE/md5sum.txt -o data/$RELEASE/md5sum.txt -z data/$RELEASE/md5sum.txt

# Consider the filenames in the md5sum file and the release notes
FILES=(`tr -s ' ' < data/$RELEASE/md5sum.txt | cut -d ' ' -f 2` release-notes.md)

if [ -d "data/$PREVIOUS" ]
then
  # Find unchanged files
  echo "Checking for unchanged files..."
  cd data/$PREVIOUS
  UNCHANGED=(`md5sum -c ../$RELEASE/md5sum.txt 2>/dev/null | grep OK |cut -d ':' -f 1  || true`)
  echo $UNCHANGED
  cd ../../

  # Hard link unchanged files
  for oldfile in "${UNCHANGED[@]}"
  do
    if [ ! -e "data/$RELEASE/$oldfile" ]
    then
      echo "Hard linking $oldfile"
      ln data/$PREVIOUS/$oldfile data/$RELEASE/$oldfile
    fi
  done
fi

# Download the items in FILES if not already present
for file in "${FILES[@]}"
do
  if [ ! -e "data/$RELEASE/$file" ]
  then
    echo "Downloading $file"
    curl $URL/$RELEASE/$file -o data/$RELEASE/$file
  fi
done


# Check the md5s for everything we downloaded except CHANGELOG.md
cd data/$RELEASE
echo "Checking MD5 hashes..."
md5sum -c md5sum.txt
cd ../../

# Make symlinks in data/ to the files in the just downloaded release folder.
for file in "${FILES[@]}"
do
  ln -sfn $RELEASE/$file data/$file
done

# make data directory unwritable in CI
if [ "$RELEASE" == "testing" ]; then
  chmod u-w data
fi
