-- Verify config_profile_item.status

BEGIN;

ALTER TABLE config.profile_item DROP COLUMN status;

ROLLBACK;
