-- Revert config_tag

BEGIN;
DROP TABLE config.tag;
COMMIT;
