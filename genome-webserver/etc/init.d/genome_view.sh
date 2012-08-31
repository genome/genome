#!/bin/sh 
### BEGIN INIT INFO
# Provides:          genome_view
# Required-Start:    $local_fs $remote_fs $network $syslog
# Required-Stop:     $local_fs $remote_fs $network $syslog
# Default-Start:     2 3 4 5
# Default-Stop:      0 1 6
# Short-Description: Start/stop genome_view fastcgi server
### END INIT INFO

# Start a Plack daemon.

PATH=/bin:/usr/bin:/sbin:/usr/sbin
USER=www-data
NAME="genome_view"
PIDFILE="/var/run/kom_fastcgi/$NAME.pid"
DAEMON="/gsc/scripts/sbin/genome_view"

# Defaults
RUN="no"
OPTIONS=""

[ -f /lib/lsb/init-functions ] && . /lib/lsb/init-functions

start()
{
    log_daemon_msg "Starting genome_view server" "$NAME"
    start-stop-daemon --start --quiet -m -b --chuid $USER --pidfile "$PIDFILE" --exec $DAEMON -- $OPTIONS

    if [ $? != 0 ]; then
        log_end_msg 1
        exit 1
    else
        log_end_msg 0
    fi
}

signal()
{

    if [ "$1" = "stop" ]; then
    SIGNAL="TERM"
        log_daemon_msg "Stopping genome_view server" "$NAME"
    else
    if [ "$1" = "reload" ]; then
        SIGNAL="HUP"
            log_daemon_msg "Reloading genome_view server" "$NAME"
    else
        echo "ERR: wrong parameter given to signal()"
        exit 1
    fi
    fi
    if [ -f "$PIDFILE" ]; then
        start-stop-daemon --stop --signal $SIGNAL --quiet --pidfile "$PIDFILE"
        if [ $? = 0 ]; then
            log_end_msg 0
        else
        SIGNAL="KILL"
        start-stop-daemon --stop --signal $SIGNAL --quiet --pidfile "$PIDFILE"
            if [ $? != 0 ]; then
                log_end_msg 1
                [ $2 != 0 ] || exit 0
            else
            rm "$PIDFILE"
                log_end_msg 0
            fi
        fi
    if [ "$SIGNAL" = "KILL" ]; then
        rm -f "$PIDFILE"
        fi
    else
        log_end_msg 0
    fi
}

case "$1" in
    start)
    start
    ;;

    force-start)
    start
    ;;

    stop)
        signal stop 0
    ;;

    force-stop)
    signal stop 0
    ;;

    reload)
    signal reload 0
    ;;

    force-reload|restart)
    signal stop 1
    sleep 2
    start
    ;;

    *)
    echo "Usage: /etc/init.d/$NAME {start|force-start|stop|force-stop|reload|restart|force-reload}"
    exit 1
    ;;
esac

exit 0

