-- Deploy instrument_data
-- requires: instrument_schema
-- requires: instrument_fragment_library

BEGIN;

CREATE TABLE IF NOT EXISTS instrument.data (
    id character varying(64) NOT NULL,
    subclass_name character varying(64) NOT NULL,
    sequencing_platform character varying(64) NOT NULL,
    library_id character varying(64) NOT NULL,
    source_name character varying(64),
    subset_name character varying(64),
    run_name character varying(64),
    CONSTRAINT data_pkey PRIMARY KEY (id),
    CONSTRAINT data_library_id_fkey FOREIGN KEY (library_id) REFERENCES instrument.fragment_library(library_id)
);

COMMIT;
