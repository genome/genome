-- Verify config_tag

BEGIN;

SELECT *
FROM config.tag
WHERE FALSE;

ROLLBACK;
