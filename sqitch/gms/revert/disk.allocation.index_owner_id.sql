-- Revert disk.allocation.index_owner_id

BEGIN;

DROP INDEX disk.d_a_owner_id;

COMMIT;
