-- Deploy model_model
-- requires: model_schema
-- requires: model_processing_profile

BEGIN;

CREATE TABLE IF NOT EXISTS model.model (
    genome_model_id character varying(32) NOT NULL,
    id character varying(32),
    name character varying(255) NOT NULL,
    sample_name character varying(255),
    processing_profile_id character varying(32) NOT NULL,
    data_directory character varying(1000),
    comparable_normal_model_id character varying(32),
    current_running_build_id character varying(32),
    last_complete_build_id character varying(32),
    user_name character varying(64),
    creation_date timestamp(6) without time zone,
    auto_assign_inst_data boolean,
    auto_build_alignments boolean,
    subject_id character varying(32) NOT NULL,
    subject_class_name character varying(500) NOT NULL,
    keep boolean DEFAULT false,
    build_granularity_unit character varying(20),
    build_granularity_value integer,
    limit_inputs_to_id character varying(1000),
    is_default boolean,
    subclass_name character varying(255),
    auto_build boolean,
    build_requested boolean,
    CONSTRAINT model_pkey PRIMARY KEY (genome_model_id),
    CONSTRAINT gm_pksni UNIQUE (genome_model_id, subclass_name),
    CONSTRAINT uniq_m_m_idun UNIQUE (genome_model_id, user_name),
    CONSTRAINT gm_gmi_sn_ppi UNIQUE (genome_model_id, subclass_name, processing_profile_id),
    CONSTRAINT model_processing_profile_id_fkey FOREIGN KEY (processing_profile_id) REFERENCES model.processing_profile(id)
);

COMMIT;
