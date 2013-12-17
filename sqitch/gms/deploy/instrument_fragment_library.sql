-- Deploy instrument_fragment_library
-- requires: instrument_schema
-- requires: subject_subject

BEGIN;

CREATE TABLE IF NOT EXISTS instrument.fragment_library (
    library_id character varying(64) NOT NULL,
    full_name character varying(64) NOT NULL,
    sample_id character varying(64) NOT NULL,
    library_insert_size character varying(64),
    original_insert_size character varying(64),
    protocol character varying(64),
    transcript_strand character varying(16),
    CONSTRAINT transcript_strand_check CHECK (((transcript_strand)::text = ANY ((ARRAY[NULL::character varying, 'unstranded'::character varying, 'firststrand'::character varying, 'secondstrand'::character varying])::text[]))),
    CONSTRAINT fragment_library_pkey PRIMARY KEY (library_id),
    CONSTRAINT fragment_library_sample_id_fkey FOREIGN KEY (sample_id) REFERENCES subject.subject(subject_id)
);

COMMIT;
