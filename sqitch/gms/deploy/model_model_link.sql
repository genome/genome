-- Deploy model_model_link
-- requires: model_model

BEGIN;

CREATE TABLE IF NOT EXISTS model.model_link (
    to_model_id character varying(32) NOT NULL,
    from_model_id character varying(32) NOT NULL,
    role character varying(56) DEFAULT 'member'::character varying NOT NULL,
    CONSTRAINT model_link_pkey PRIMARY KEY (to_model_id, from_model_id),
    CONSTRAINT model_link_from_model_id_fkey FOREIGN KEY (from_model_id) REFERENCES model.model(genome_model_id),
    CONSTRAINT model_link_to_model_id_fkey FOREIGN KEY (to_model_id) REFERENCES model.model(genome_model_id)
);

COMMIT;
