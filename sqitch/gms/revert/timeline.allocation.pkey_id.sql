-- Revert timeline.allocation.pkey_id

BEGIN;
    ALTER TABLE timeline.allocation DROP CONSTRAINT allocation_pkey;
COMMIT;
