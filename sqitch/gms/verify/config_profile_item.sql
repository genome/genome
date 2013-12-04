-- Verify config_profile_item

BEGIN;

SELECT *
FROM config.profile_item
WHERE FALSE;

ROLLBACK;
