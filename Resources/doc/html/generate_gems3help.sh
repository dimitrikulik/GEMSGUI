#!/bin/bash
# Re-generate gems3help.qhc and gems3help.qch from Qt Help project configuration

QHELPGEN=""
if [ -n "$1" ] && [ -x "$1" ]; then
    QHELPGEN="$1"
elif command -v qhelpgenerator &>/dev/null; then
    QHELPGEN="qhelpgenerator"
elif [ -n "$QTDIR" ] && [ -x "$QTDIR/libexec/qhelpgenerator" ]; then
    QHELPGEN="$QTDIR/libexec/qhelpgenerator"
elif [ -n "$QTDIR" ] && [ -x "$QTDIR/bin/qhelpgenerator" ]; then
    QHELPGEN="$QTDIR/bin/qhelpgenerator"
elif [ -x "/Users/kulik/Qt/6.11.2/macos/libexec/qhelpgenerator" ]; then
    QHELPGEN="/Users/kulik/Qt/6.11.2/macos/libexec/qhelpgenerator"
elif [ -x "/home/sveta/Qt/6.8.2/gcc_64/libexec/qhelpgenerator" ]; then
    QHELPGEN="/home/sveta/Qt/6.8.2/gcc_64/libexec/qhelpgenerator"
fi

if [ -z "$QHELPGEN" ]; then
    echo "Error: qhelpgenerator not found. Pass path to qhelpgenerator as argument or set QTDIR."
    exit 1
fi

cd "$(dirname "$0")" || exit 1
echo "Generating Qt Help with: $QHELPGEN"
"$QHELPGEN" gems3helpconfig.qhcp -o gems3help.qhc
