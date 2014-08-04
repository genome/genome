-- Revert config_tag_profile_item

BEGIN;
DROP TABLE config.tag_profile_item;
COMMIT;
