#! /bin/bash

echo "Some STDERR Text" 1>&2
echo "Some STDOUT Text"

sleep 1

echo "STDERR: Launching Process /this/is/a/test/ (foo)" 1>&2

sleep 1

echo "Some more STDERR Text" 1>&2
echo "Some more STDOUT Text"
