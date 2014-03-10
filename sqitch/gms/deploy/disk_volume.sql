-- Deploy disk_volume
-- requires: disk_schema

BEGIN;

CREATE TABLE IF NOT EXISTS disk.volume (
    id character varying(255) NOT NULL,
    hostname character varying(255) NOT NULL,
    physical_path character varying(255) NOT NULL,
    mount_path character varying(255) NOT NULL,
    total_kb bigint NOT NULL,
    unallocated_kb bigint NOT NULL,
    disk_status character varying(15) NOT NULL,
    can_allocate boolean DEFAULT true NOT NULL,
    doubles_space boolean DEFAULT false NOT NULL,
    CONSTRAINT volume_pkey PRIMARY KEY (id)
);


COMMIT;
