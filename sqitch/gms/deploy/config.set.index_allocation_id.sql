-- Deploy config.set.index_allocation_id
-- requires: config_set

BEGIN;

CREATE INDEX c_s_allocation_id_index ON config.set USING btree (allocation_id);

COMMIT;
