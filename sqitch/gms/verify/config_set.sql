-- Verify config_set

BEGIN;

SELECT id, created_at, updated_at, allocation_id
FROM config.set
WHERE FALSE;

ROLLBACK;
