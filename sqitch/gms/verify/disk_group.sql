-- Verify disk_group

BEGIN;

SELECT id, name, permissions, sticky, subdirectory, unix_uid, unix_gid
FROM disk."group"
WHERE FALSE;

ROLLBACK;
