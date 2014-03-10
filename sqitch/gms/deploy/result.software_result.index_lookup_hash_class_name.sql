-- Deploy result.software_result.lookup_hash_class_name
-- requires: result_software_result

BEGIN;

CREATE INDEX lookup_hash_class_name_idx on result.software_result using btree (lookup_hash, class_name);

COMMIT;
