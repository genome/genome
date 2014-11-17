-- Deploy process_input
-- requires: process_process

BEGIN;

CREATE TABLE IF NOT EXISTS process.input (
    id character varying(64) PRIMARY KEY,
    process_id character varying(64) NOT NULL,
    value_class_name text,
    value_id text NOT NULL,
    name text NOT NULL,
    array_index int NOT NULL,
    CONSTRAINT input_process_fkey FOREIGN KEY (process_id) REFERENCES process.process(id),
    CONSTRAINT process_input_unique_name_index UNIQUE (process_id, name, array_index)
);

REVOKE ALL ON TABLE process.input FROM PUBLIC;
GRANT ALL ON TABLE process.input TO genome;
GRANT SELECT,INSERT ON TABLE process.input TO "gms-user";

COMMIT;
