-- Deploy disk_volume_group_bridge
-- requires: disk_volume
-- requires: disk_group

BEGIN;

CREATE TABLE IF NOT EXISTS disk.volume_group_bridge (
    volume_id character varying(255) NOT NULL,
    group_id character varying(255) NOT NULL,
    CONSTRAINT volume_group_bridge_pkey PRIMARY KEY (volume_id, group_id)
);

COMMIT;
