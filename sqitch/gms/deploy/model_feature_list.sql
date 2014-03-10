-- Deploy model_feature_list
-- requires: model_schema
-- requires: model_build

BEGIN;

CREATE TABLE IF NOT EXISTS model.feature_list (
    id character varying(64) NOT NULL,
    name character varying(200) NOT NULL,
    source character varying(64),
    format character varying(64) NOT NULL,
    reference_id character varying(32),
    subject_id character varying(32),
    file_id bigint,
    file_content_hash character varying(32) NOT NULL,
    content_type character varying(255),
    description character varying(255),
    CONSTRAINT feature_list_pkey PRIMARY KEY (id),
    CONSTRAINT feature_list_name_key UNIQUE (name),
    CONSTRAINT feature_list_reference_id_fkey FOREIGN KEY (reference_id) REFERENCES model.build(build_id),
    CONSTRAINT feature_list_subject_id_fkey FOREIGN KEY (subject_id) REFERENCES model.build(build_id)
);

COMMIT;
