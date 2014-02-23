-- Deploy config.subject_mapping
-- requires: config_analysis_project
-- requires: subject_pairing

BEGIN;

CREATE TABLE IF NOT EXISTS config.subject_mapping (
    id character varying(64) NOT NULL PRIMARY KEY,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    created_by character varying(255) NOT NULL,
    analysis_project_id character varying(64) NOT NULL REFERENCES config.analysis_project(id)
);

CREATE TABLE IF NOT EXISTS config.subject_mapping_input (
    id character varying(64) NOT NULL PRIMARY KEY,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    created_by character varying(255) NOT NULL,
    subject_mapping_id character varying(64) NOT NULL REFERENCES config.subject_mapping(id),
    key text NOT NULL,
    value text NOT NULL
);

CREATE TABLE IF NOT EXISTS config.subject_mapping_subject (
    id character varying(64) NOT NULL PRIMARY KEY,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    created_by character varying(255) NOT NULL,
    subject_mapping_id character varying(64) NOT NULL REFERENCES config.subject_mapping(id),
    label text NOT NULL,
    subject_id character varying(64) REFERENCES subject.subject(subject_id)
);

REVOKE ALL ON TABLE config.subject_mapping FROM PUBLIC;
REVOKE ALL ON TABLE config.subject_mapping FROM genome;
GRANT ALL ON TABLE config.subject_mapping TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.subject_mapping TO "gms-user";

REVOKE ALL ON TABLE config.subject_mapping_input FROM PUBLIC;
REVOKE ALL ON TABLE config.subject_mapping_input FROM genome;
GRANT ALL ON TABLE config.subject_mapping_input TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.subject_mapping_input TO "gms-user";

REVOKE ALL ON TABLE config.subject_mapping_subject FROM PUBLIC;
REVOKE ALL ON TABLE config.subject_mapping_subject FROM genome;
GRANT ALL ON TABLE config.subject_mapping_subject TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.subject_mapping_subject TO "gms-user";

DROP TABLE IF EXISTS subject.pairing;

COMMIT;
