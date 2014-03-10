-- Deploy disk.alloction.index_allocation_path
-- requires: disk_allocation

BEGIN;

CREATE INDEX idx_allocation_path_like ON disk.allocation USING btree (allocation_path varchar_pattern_ops);

COMMIT;
