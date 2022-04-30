#!/bin/bash
versionfile="../version.txt"
citationfile="../CITATION.cff"
do_release() {
    git status;
    new_version=$(awk -F '.' '{printf("%i.%i.%i\n", $1, $2, $3+1)}' < $versionfile)
    read -p "New version (^C to cancel, empty for $new_version: " version_given
    if [ ! -z "$version_given" ]; then
        new_version="$version_given";
    fi
    echo "You have chosen: $new_version"
    echo "$new_version" > "$versionfile"
    ./citation.sh
    git add "$versionfile" "$citationfile"
    git commit -m "Version bump to $new_version"
    git tag "v$new_version"
    echo "You can now push (I won't do it)"
}

if [ ! -z "$(git status -s)" ]; then
    echo "There are modified or untracked files present, clean up first!"
    git status -s;
    echo "Release aborted!"
    exit 1;
fi

echo "Change version number, commit and tag the repository."
echo -n "Current version is: "
cat "$versionfile"
read -p "Do you wish to continue? " reply
case $reply in 
    [Yy]* ) do_release;;
esac
