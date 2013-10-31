-- Verify web_talk

BEGIN;

SELECT id, user_id, command_class, stdout_pathname, stderr_pathname,
    status, time_submitted, time_started, time_finished
FROM web.task
WHERE FALSE;

ROLLBACK;
