-- Deploy config_profile_item.status
-- requires: config_profile_item

BEGIN;

ALTER TABLE config.profile_item ADD COLUMN status TEXT;

COMMIT;
