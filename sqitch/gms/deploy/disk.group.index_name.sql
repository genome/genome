-- Deploy disk.group.name
-- requires: disk_group

BEGIN;

CREATE INDEX group_name_index on disk."group" using btree (name);

COMMIT;
