-- Deploy disk.allocation.index_owner_class_name_owner_id
-- requires: disk_allocation

BEGIN;

CREATE INDEX allocation_owner_class_id_index ON disk.allocation USING btree (owner_class_name, owner_id);

COMMIT;
