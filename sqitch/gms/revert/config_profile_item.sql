-- Revert config_profile_item

BEGIN;
DROP TABLE config.profile_item;
COMMIT;
