-- Deploy disk.file_summary.index_allocation_id
-- requires: disk_file_summary

BEGIN;

CREATE INDEX file_summary_allocation_id_idx ON disk.file_summary USING btree (allocation_id);

COMMIT;
