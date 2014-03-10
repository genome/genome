-- Deploy model_build_link
-- requires: model_build

BEGIN;

CREATE TABLE IF NOT EXISTS model.build_link (
    to_build_id character varying(64) NOT NULL,
    from_build_id character varying(64) NOT NULL,
    role character varying(56) NOT NULL,
    CONSTRAINT build_link_pkey PRIMARY KEY (to_build_id, from_build_id),
    CONSTRAINT build_link_from_build_id_fkey FOREIGN KEY (from_build_id) REFERENCES model.build(build_id),
    CONSTRAINT build_link_to_build_id_fkey FOREIGN KEY (to_build_id) REFERENCES model.build(build_id)
);

COMMIT;
