-- Deploy disk.volume.mount_path
-- requires: disk_volume

BEGIN;

CREATE INDEX volume_mount_path_index on disk.volume using btree (mount_path);

COMMIT;
