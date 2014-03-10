-- Revert config_set

BEGIN;

DROP TABLE IF EXISTS config.set;

COMMIT;
