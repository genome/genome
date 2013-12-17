-- Revert config_analysismenu_item

BEGIN;
DROP TABLE config.analysismenu_item;
COMMIT;
