-- Deploy disk.allocation_permission_scheme
-- requires: disk_allocation

BEGIN;

CREATE TABLE disk.allocation_permission_scheme (
    id character varying(64) NOT NULL,
    uid integer NOT NULL,
    gid integer NOT NULL,
    min_mode integer NOT NULL,
    max_mode integer NOT NULL,
    file_min_mask integer NOT NULL,
    CONSTRAINT allocation_permission_scheme_ok PRIMARY KEY (id)
);

ALTER TABLE disk.allocation ADD COLUMN permission_scheme_id character varying(64);
ALTER TABLE disk.allocation ADD CONSTRAINT allocation_permission_scheme_fk
    FOREIGN KEY (permission_scheme_id)
    REFERENCES disk.allocation_permission_scheme(id);

COMMIT;
