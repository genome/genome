-- Deploy disk_allocation
-- requires: disk_schema

BEGIN;

CREATE TABLE IF NOT EXISTS disk.allocation (
    id character varying(255) NOT NULL,
    allocation_path character varying(4000) NOT NULL,
    disk_group_name character varying(40) NOT NULL,
    group_subdirectory character varying(255) NOT NULL,
    mount_path character varying(255) NOT NULL,
    kilobytes_requested numeric(20,0) NOT NULL,
    kilobytes_used numeric(20,0),
    owner_class_name character varying(255),
    owner_id character varying(255),
    creation_time timestamp(6) without time zone,
    reallocation_time timestamp(6) without time zone,
    original_kilobytes_requested numeric(20,0),
    archivable boolean,
    preserved boolean,
    kilobytes_used_time timestamp without time zone,
    archive_after_time timestamp without time zone,
    status character varying(40),
    CONSTRAINT allocation_pkey PRIMARY KEY (id),
    CONSTRAINT d_a_uniq_mount_path UNIQUE (mount_path, group_subdirectory, allocation_path),
    CONSTRAINT uniq_d_a_apid UNIQUE (allocation_path, id)
);

COMMIT;
