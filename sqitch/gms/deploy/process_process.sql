-- Deploy process_process
-- requires: process_schema

BEGIN;

CREATE TABLE IF NOT EXISTS process.process (
    id character varying(64) PRIMARY KEY,
    disk_allocation_id character varying(255) NOT NULL UNIQUE,
    subclass_name text NOT NULL,
    created_by text NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    started_at timestamp(6) without time zone,
    ended_at timestamp(6) without time zone,
    status text,
    software_revision text,
    CONSTRAINT process_disk_allocation_fkey FOREIGN KEY (disk_allocation_id) REFERENCES disk.allocation(id)
);

REVOKE ALL ON TABLE process.process FROM PUBLIC;
GRANT ALL ON TABLE process.process TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE process.process TO "gms-user";

COMMIT;
