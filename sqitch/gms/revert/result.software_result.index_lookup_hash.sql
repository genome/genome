-- Revert result.software_result.index_lookup_hash

BEGIN;

DROP INDEX result.lookup_hash_idx;

COMMIT;
