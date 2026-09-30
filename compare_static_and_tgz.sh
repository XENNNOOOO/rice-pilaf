tar -tf static.tgz | sed 's#/$##' > tgz
find static > fs
sort -o fs fs & sort -o tgz tgz & wait
diff -uw fs tgz > diff.txt