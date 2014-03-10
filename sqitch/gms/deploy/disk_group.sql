-- Deploy disk_group
-- requires: disk_schema

BEGIN;

CREATE TABLE IF NOT EXISTS disk."group" (
    id character varying(255) NOT NULL,
    name character varying(40) NOT NULL,
    permissions integer NOT NULL,
    sticky boolean NOT NULL,
    subdirectory character varying(255) NOT NULL,
    unix_uid integer NOT NULL,
    unix_gid integer NOT NULL,
    CONSTRAINT group_pkey PRIMARY KEY (id)
);

COMMIT;
