-- Deploy disk.allocation.index_owner_id
-- requires: disk_allocation

BEGIN;

CREATE INDEX d_a_owner_id ON disk.allocation USING btree (owner_id);

COMMIT;
