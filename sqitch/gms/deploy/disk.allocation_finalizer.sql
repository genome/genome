-- Deploy disk.allocation_finalizer
-- requires: disk_allocation

BEGIN;

CREATE TABLE disk.allocation_finalizer (
    id character varying(64) NOT NULL,
    uid integer NOT NULL,
    gid integer NOT NULL,
    min_mode integer NOT NULL,
    max_mode integer NOT NULL,
    file_min_mask integer NOT NULL,
    CONSTRAINT allocation_finalizer_ok PRIMARY KEY (id)
);

ALTER TABLE disk.allocation ADD COLUMN finalizer_id character varying(64);
ALTER TABLE disk.allocation ADD CONSTRAINT allocation_finalizer_fk
    FOREIGN KEY (finalizer_id)
    REFERENCES disk.allocation_finalizer(id);

COMMIT;
