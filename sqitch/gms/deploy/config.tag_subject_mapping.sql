-- Deploy config.tag_subject_mapping
-- requires: config_tag
-- requires: config.subject_mapping

BEGIN;

CREATE TABLE IF NOT EXISTS config.tag_subject_mapping (
    id character varying(64) PRIMARY KEY,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    created_by character varying(255) NOT NULL,
    tag_id character varying(64) NOT NULL REFERENCES config.tag(id),
    subject_mapping_id character varying(64) NOT NULL REFERENCES config.subject_mapping(id)
);

REVOKE ALL ON TABLE config.tag_subject_mapping FROM PUBLIC;
GRANT ALL ON TABLE config.tag_subject_mapping TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.tag_subject_mapping TO "gms-user";

COMMIT;
