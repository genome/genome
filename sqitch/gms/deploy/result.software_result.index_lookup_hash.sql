-- Deploy result.software_result.lookup_hash
-- requires: result_software_result

BEGIN;

CREATE INDEX lookup_hash_idx on result.software_result using btree (lookup_hash);

COMMIT;
