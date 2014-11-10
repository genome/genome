-- Deploy process_status_event
-- requires: process_schema

BEGIN;

CREATE TABLE IF NOT EXISTS process.status_event (
    id character varying(64) PRIMARY KEY,
    process_id character varying(64) NOT NULL,
    old_status text,
    new_status text NOT NULL,
    timestamp timestamp(6) without time zone NOT NULL,
    CONSTRAINT status_event_process_fkey FOREIGN KEY (process_id) REFERENCES process.process(id)
);

REVOKE ALL ON TABLE process.status_event FROM PUBLIC;
GRANT ALL ON TABLE process.status_event TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE process.status_event TO "gms-user";

COMMIT;
