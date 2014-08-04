-- Verify config_tag_profile_item

BEGIN;

SELECT *
FROM config.tag_profile_item
WHERE FALSE;

ROLLBACK;
